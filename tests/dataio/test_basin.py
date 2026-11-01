"""Pins for mls_emergence.dataio.basin (probe spec rows 3 and 7, pilot step 2).

Expected values never come through the function under test: the counts and
coordinates are compared with the Codex audit's independently saved copies
(codex/basin_counts.csv, codex/basin_coords.csv, written by codex/audit_checks.py
from the project loaders on 2026-10-16), the distance is compared with a
great-circle formula written here, and the membership rule is re-applied from
its ingredients rather than read back from the file it produced.
"""
from __future__ import annotations

import csv
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

from mls_emergence.dataio import basin

ROOT = Path(__file__).resolve().parents[2]
AUDIT_COUNTS = ROOT / "codex" / "basin_counts.csv"
AUDIT_COORDS = ROOT / "codex" / "basin_coords.csv"


@pytest.fixture(scope="module")
def loaded():
    counts, coords = basin.load_curated()
    return counts, coords


def _read_audit(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]
    return header, body


def test_loader_reproduces_the_audits_counts_exactly(loaded):
    counts, _ = loaded
    header, body = _read_audit(AUDIT_COUNTS)
    assert header[1:] == list(counts.columns) == basin.DECORATED_TYPES
    assert [r[0] for r in body] == list(counts.index)
    audit = np.array([[int(v) for v in r[1:]] for r in body])
    assert counts.shape == (28, 10)
    assert np.array_equal(counts.to_numpy(float), audit.astype(float))
    assert int(audit.sum()) == 14101


def test_loader_reproduces_the_audits_coordinates(loaded):
    _, coords = loaded
    header, body = _read_audit(AUDIT_COORDS)
    lat_i, lon_i = header.index("Latitude"), header.index("Longitude")
    assert [r[0] for r in body] == list(coords.index)
    audit = np.array([[float(r[lat_i]), float(r[lon_i])] for r in body])
    assert np.allclose(coords[["Latitude", "Longitude"]].to_numpy(float), audit, atol=1e-9, rtol=0)


def test_membership_rule_reapplied_from_its_ingredients_gives_the_committed_list():
    """The committed list equals the phase rule applied afresh (rule 19).

    Ingredients: the whole-LMV curated set (script 12), Phillips's phase
    labels (script 36) and the three constants of script 16. The rule is
    re-stated here rather than called, so that a change to either side shows.
    """
    sys.path.insert(0, str(ROOT / "analyses"))
    m16 = importlib.import_module("16_basin_membership")
    m36 = importlib.import_module("36_canonical_phase_map")
    m12 = importlib.import_module("12_sensitivity_grid")
    counts, coords = m12.load_curated_full()
    coords = coords.dropna()
    names = [str(i) for i in coords.index]
    labels, _derived = m36.assign_primary_phases(names, coords[["Latitude", "Longitude"]].to_numpy(float))
    labels = np.asarray(labels)
    sherds = counts.reindex(coords.index).sum(axis=1).to_numpy(float)
    keep = np.isin(labels, m16.ST_FRANCIS_PHASES) & (sherds >= m16.MIN_DECORATED)
    for phase in m16.ST_FRANCIS_PHASES:
        if 0 < int((keep & (labels == phase)).sum()) < m16.MIN_PHASE_MEMBERS:
            keep &= labels != phase
    regenerated = sorted(np.asarray(names)[keep])
    assert regenerated == sorted(basin.basin_members("curated"))
    assert len(regenerated) == 28
    # The rule discriminates: without the sherd threshold the set is larger.
    assert int((np.isin(labels, m16.ST_FRANCIS_PHASES)).sum()) > 28


def _great_circle_km(c):
    lat, lon = np.deg2rad(c[:, 0]), np.deg2rad(c[:, 1])
    h = (np.sin((lat[:, None] - lat[None, :]) / 2) ** 2
         + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin((lon[:, None] - lon[None, :]) / 2) ** 2)
    return 6371.0088 * 2 * np.arcsin(np.sqrt(h))


def test_geo_km_is_within_one_percent_of_great_circle_and_is_not_great_circle(loaded):
    _, coords = loaded
    c = coords.to_numpy(float)
    D = basin.geo_km(c)
    G = _great_circle_km(c)
    off = ~np.eye(len(c), dtype=bool)
    rel = np.abs(D[off] - G[off]) / G[off]
    assert rel.max() < 0.01, f"equirectangular error {rel.max():.4f} exceeds one percent"
    # The approximation is not exact: a great-circle implementation would make
    # this zero, so the assertion records that geo_km is the equirectangular form.
    assert rel.max() > 1e-6
    assert np.allclose(D, D.T) and np.all(np.diag(D) == 0)


def test_geo_km_equals_distances_between_km_xy_rows(loaded):
    _, coords = loaded
    c = coords.to_numpy(float)
    p = basin.km_xy(c)
    D2 = np.linalg.norm(p[:, None] - p[None, :], axis=-1)
    assert np.max(np.abs(basin.geo_km(c) - D2)) < 1e-9


@pytest.mark.parametrize("bad", [
    np.array([35.0, -90.0]),                       # 1-D
    np.array([[35.0, -90.0, 1.0]]),                # three columns
    np.array([[35.0, np.nan], [35.1, -90.0]]),     # NaN
    np.array([[135.0, -90.0], [35.1, -90.0]]),     # latitude out of range
])
def test_geo_km_and_km_xy_refuse_out_of_domain_input(bad):
    with pytest.raises(ValueError):
        basin.geo_km(bad)
    with pytest.raises(ValueError):
        basin.km_xy(bad)


def test_fingerprint_is_stable_and_reports_the_measured_geometry(loaded):
    counts, coords = loaded
    a = basin.data_fingerprint(counts, coords)
    b = basin.data_fingerprint(counts, coords)
    assert a == b
    assert a["n_assemblages"] == 28 and a["n_classes"] == 10 and a["n_sherds"] == 14101
    # Geometry against the great-circle formula, not against a quoted number.
    G = _great_circle_km(coords.to_numpy(float))
    assert abs(a["max_pairwise_km"] - G.max()) / G.max() < 0.01
    assert a["provisional_coordinates"] == ["Cummins"]
    assert basin.provisional_coordinates(["Parkin", "Cummins", "Nowhere"]) == ["Cummins"]
    # A changed count changes the hash (the fingerprint is not a constant).
    c2 = counts.copy(); c2.iloc[0, 0] += 1
    assert basin.data_fingerprint(c2, coords)["sha256_counts"] != a["sha256_counts"]
    with pytest.raises(ValueError):
        basin.data_fingerprint(counts.iloc[::-1], coords)
