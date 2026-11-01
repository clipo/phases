"""Placements from geography (A3 design, section 3)."""
import csv
from pathlib import Path

import numpy as np

from mls_emergence.dataio.basin import load_curated
from mls_emergence.probe import placement as pl

ROOT = Path(__file__).resolve().parents[2]


def test_line_regions_split_known_points():
    pts = np.array([[35.0, -90.5], [35.0, -90.0], [35.5, -90.5], [35.5, -90.0]])   # SW, SE, NW, NE
    # Vertical line (angle 0: projection is easting) at offset 0 (the centroid): west sites region 1.
    assert pl.line_regions(pts, 0, 0.0).tolist() == [1, 2, 1, 2]
    # Horizontal line (angle 90: projection is northing): south sites region 1.
    assert pl.line_regions(pts, 90, 0.0).tolist() == [1, 1, 2, 2]


def test_candidate_set_has_no_duplicates_and_respects_minimum_side():
    _, coords = load_curated(); xy = coords.to_numpy(float)
    cands = pl.candidate_lines(xy)
    keys = [tuple(c["regions"].tolist()) for c in cands]
    assert len(keys) == len(set(keys))
    assert all(min(c["sides"]) >= pl.MIN_SIDE for c in cands)
    assert 10 <= len(cands) <= 30


def test_kmeans2_regions_are_the_screens_two_clusters_and_split_walls_off():
    counts, coords = load_curated(); xy = coords.to_numpy(float)
    reg = pl.kmeans2_regions(xy)
    assert set(reg.tolist()) == {1, 2}
    with open(ROOT / "codex" / "basin_membership_check.csv", newline="") as f:
        area = {r["assemblage"]: r["phase"] for r in csv.DictReader(f)}
    names = list(counts.index)
    # The eastern cluster (label 2, larger easting) must contain every Walls site.
    walls = [i for i, n in enumerate(names) if area[n] == "Walls"]
    assert all(reg[i] == 2 for i in walls)


def test_best_line_for_recovers_a_line_partition():
    _, coords = load_curated(); xy = coords.to_numpy(float)
    cands = pl.candidate_lines(xy)
    target = cands[7]["regions"]
    b = pl.best_line_for(target, cands)
    assert b["agreement"] == 1.0 and b["angle"] == cands[7]["angle"]
