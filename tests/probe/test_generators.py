"""Gabriel edges and graph copying weights (D21)."""
import csv
from pathlib import Path

import numpy as np
import pytest

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.probe import generators as g

ROOT = Path(__file__).resolve().parents[2]


def test_gabriel_edges_on_a_hand_checkable_configuration():
    # Four points on a square plus its center: the center is a Gabriel neighbor
    # of every corner; opposite corners are not neighbors (the center sits on
    # the circle's diameter); adjacent corners are neighbors.
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]], float)
    D = np.linalg.norm(pts[:, None] - pts[None, :], axis=-1)
    E = {tuple(e) for e in g.gabriel_edges(D)}
    assert E == {(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)}
    with pytest.raises(ValueError):
        g.gabriel_edges(D[:3, :2])


def test_gabriel_edges_on_the_28_sites_reproduce_the_design_measurement():
    """41 edges, 8 crossing Phillips's areas, as the committed design script found (D14)."""
    counts, coords = load_curated()
    D = geo_km(coords.to_numpy(float))
    E = g.gabriel_edges(D)
    assert len(E) == 41
    with open(ROOT / "codex" / "basin_membership_check.csv", newline="") as f:
        area = {r["assemblage"]: r["phase"] for r in csv.DictReader(f)}
    names = list(counts.index)
    assert sum(area[names[i]] != area[names[j]] for i, j in E) == 8


def test_graph_weights_are_row_stochastic_equal_on_edges_and_respect_a_barrier():
    E = np.array([[0, 1], [1, 2], [2, 3]])
    w = g.graph_weights(E, 5)            # node 4 has no edge: copies itself
    assert np.allclose(w.sum(1), 1)
    assert w[1, 0] == 0.5 and w[1, 2] == 0.5 and w[0, 1] == 1.0 and w[4, 4] == 1.0
    labels = np.array([0, 0, 1, 1, 1])
    wb = g.graph_weights(E, 5, labels=labels, leak=0.1)
    assert np.isclose(wb[1, 2], 0.1 / 1.1) and np.isclose(wb[1, 0], 1.0 / 1.1)   # cross-label edge down-weighted
    wc = g.graph_weights(E, 5, labels=labels, leak=0.0)
    assert wc[1, 2] == 0 and wc[1, 0] == 1.0                                     # a full barrier cuts the edge
    with pytest.raises(ValueError):
        g.graph_weights(E, 5, leak=1.5)


def test_family_weights_dispatch():
    pts = np.array([[35.0, -90.0], [35.1, -90.0], [35.0, -90.2]])
    D = geo_km(pts)
    assert np.allclose(g.family_weights("exponential", 12.0, D).sum(1), 1)
    assert np.allclose(g.family_weights("gabriel", None, D).sum(1), 1)
    assert g.family_name("gabriel", None) == "gabriel" and g.family_name("power", 6.0) == "power-6km"
    assert len(g.FAMILIES) == 10
