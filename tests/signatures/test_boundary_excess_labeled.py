"""The phase-line boundary excess must not rest on a single pair of assemblages.

`23_phases_vs_spatial_drift.boundary_excess_labeled` averages the
within-minus-between similarity gap over distance bins. Until 2026-10-02 a bin
counted with one pair of each kind, and on the basin's river-distance matrix
the farthest bin held one within-phase pair: correcting one coordinate changed
which pair it was and flipped the sign of the statistic.

Rule 3: the probe is a record on which the two definitions disagree in SIGN.
Two groups of four sites interleave along 70 km and differ in composition, so
at short range within-group pairs are more alike than between-group pairs.
Three outliers (two of group A at opposite ends, one of group B beside the
second) put exactly one within-group pair and one between-group pair in the
farthest bin, and that within-group pair is the least alike in the record.
The expected values in the comments were worked by hand from the profiles.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
sd = importlib.import_module("23_phases_vs_spatial_drift")


def _record():
    A, B = (40.0, 40.0, 20.0), (30.0, 30.0, 40.0)        # BR similarity A-B: 160
    P, Q, R = (80.0, 0.0, 20.0), (0.0, 80.0, 20.0), (70.0, 10.0, 20.0)
    #        x (km)  group  profile
    sites = [(10, 0, A), (30, 0, A), (50, 0, A), (70, 0, A),
             (0, 1, B), (20, 1, B), (40, 1, B), (60, 1, B),
             (-600, 0, P),      # group A, far west
             (600, 0, Q),       # group A, far east: P-Q similarity 40
             (610, 1, R)]       # group B, beside Q: P-R similarity 180
    x = np.array([t[0] for t in sites], float)
    labels = np.array([t[1] for t in sites])
    counts = np.array([t[2] for t in sites], float)
    return counts, np.abs(x[:, None] - x[None, :]), labels


def test_one_far_pair_no_longer_sets_the_sign():
    counts, dist, labels = _record()
    old = sd.boundary_excess_labeled(counts, dist, labels, min_pairs=1)
    new = sd.boundary_excess_labeled(counts, dist, labels)
    # Bin gaps by hand, nearest to farthest: +45.9 (12 within, 17 between),
    # +2.2 (7, 9), +20.0 (5 within but only 3 between), -140.0 (1, 1).
    # One pair of each kind in the last bin drags the old mean below zero;
    # the new one keeps the two bins with five or more pairs of each kind.
    assert old == pytest.approx((45.882 + 2.222 + 20.0 - 140.0) / 4, abs=0.01)
    assert new == pytest.approx((45.882 + 2.222) / 2, abs=0.01)
    assert old < 0 < new


def test_bins_below_the_minimum_are_left_out_not_zeroed():
    counts, dist, labels = _record()
    # With a minimum no bin can meet, the function falls back to the unbinned
    # gap rather than returning 0 or raising.
    S = sd.similarity_matrix(counts)
    iu = np.triu_indices(len(counts), 1)
    same = labels[iu[0]] == labels[iu[1]]
    raw = float(S[iu][same].mean() - S[iu][~same].mean())
    assert sd.boundary_excess_labeled(counts, dist, labels, min_pairs=10_000) == pytest.approx(raw)


def test_minimum_must_be_positive():
    counts, dist, labels = _record()
    with pytest.raises(ValueError):
        sd.boundary_excess_labeled(counts, dist, labels, min_pairs=0)
